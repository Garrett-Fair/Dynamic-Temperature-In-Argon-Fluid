// This is main.
// This is MD simulation of many particles interacting with pair potentials

using namespace std;

#include "functions.h"
#include <chrono>
#include <random>
#include <cmath>

void compute_forces(vector<PARTICLE>&, double, double, double);
double compute_potential_energies(vector<PARTICLE>&, double, double, double);

int main(int argc, char* argv[]) 
{
  // we begin with defining the Boltzmann constant
  const double kB = 1.38e-23;	// Joules per Kelvin

  cout << "\n---Simulating fluid (Argon)---\n";

  // key descriptors of the many particle system
  double diameter = 3.405e-10;  // what is the size of each particle? (in meters)
  double mass = 6.634e-26;  // what is the mass of each particle? (in kg)
  double ljenergy = 120 * kB; // what is the strength of typical particle-particle collisions? (in Joules)

  // reduced units
  double unitlength = diameter; // unit of length is the simulation (diameter of the particle)
  double unitmass = mass; // unit of mass in the simulation (mass of the particle)
  double unitenergy = ljenergy; // unit of energy in the simulation (characteristic pair potential strength)
  double unittemperature = ljenergy/kB;

  double unittime = sqrt(unitmass * unitlength * unitlength / unitenergy); // Unit of time

  cout << "unit of length is " << unitlength << " meters" << endl;
  cout << "unit of mass is " << unitmass << " kilograms" << endl;
  cout << "unit of energy is " << unitenergy << " Joules" << endl;
  cout << "unit of time is " << unittime << " seconds" << endl;
  cout << "unit of temperature is " << unittemperature << " Kelvin" << endl;
  cout << "\n";

  double reduced_diameter = diameter / unitlength;
  double reduced_mass = mass / unitmass;
  double reduced_ljenergy = ljenergy / unitenergy;

  int Nparticles = 216;
  cout << "initializing 216 particles on a lattice (simple cubic). density will be changed by changing box size" << endl;

  double density = 0.8442;
  cout << "enter density (in reduced units; suggested 0.8442) " << endl;
  cin >> density;

  double reduced_temperature = 1;
  cout << "enter temperature (in reduced units; suggested 1) " << endl;
  cin >> reduced_temperature;
  double temperature = reduced_temperature * unittemperature; // in Kelvins



  double temperature_change = 1;
  cout << "enter temperature change (in reduced units) " << endl;
  cin >> temperature_change;
  //temperature = temperature_change * unittemperature; //in Kelvins








  double initial_temperature = 0.728;   // in reduced units

  double distance_cutoff = 2.5; // in units of the diameter (chosen as unitlength)

  double boxL = pow(Nparticles/density,1.0/3.0); // box length
  cout << "cubic box edge length (in reduced units) calculated using " << Nparticles << " particles and " << density << " density: " << boxL << endl;

  if (boxL <= 2*distance_cutoff)
  {
      cout << "the box edge length is not large enough to make PBC work" << endl;
      cout << "not proceeding further. increase number of particles to resolve the issue. choose a cubic number." << endl;
      return 0;
  }

  // Different parts of the system
  vector<PARTICLE> particle;		// all particles in the system

  initialize_particle_positions(particle, density, Nparticles, reduced_diameter, reduced_mass, boxL);

  initialize_particle_velocities(particle, initial_temperature);

  // output to screen
  cout << "Number of particles inside the box: " << particle.size() << endl;
  cout << "density of the fluid is " << particle.size()/(boxL*boxL*boxL) << " and in kg/m^3 " << density*mass/(pow(diameter,3)) << endl;
  cout << "temperature of the fluid is " << reduced_temperature << " and in Kelvin " << temperature << endl;

  cout << "\n";

  // initial energies and forces computation
  double totalke = 0.0;
  for (unsigned int i = 0; i < particle.size(); i++)
  {
      particle[i].kinetic_energy();
      totalke += particle[i].ke;
  }

  double totalpe = compute_potential_energies(particle, reduced_ljenergy, distance_cutoff, boxL);
  compute_forces(particle, reduced_ljenergy, distance_cutoff, boxL);

  // create files for storing movie and energies

  char file_movie[200], file_energy[200], file_temperature[200], file_velocity[200], file_velocity_all[200];

  // file movie

  sprintf(file_movie, "movie_rho%f_T%f.out", density, reduced_temperature);
  ofstream list_propagation(file_movie, ios::out); // create a file to store and visualize 3D data
  make_movie(0,particle,list_propagation, boxL);

  // file ke, pe, and total energies of all particles
  sprintf(file_energy, "energy_rho%f_T%f.out", density, reduced_temperature);
  ofstream output_energy(file_energy, ios::out);

  // file ke, pe, and total energies of all particles
  sprintf(file_temperature, "temperature_rho%f_T%f.out", density, reduced_temperature);
  ofstream output_temperature(file_temperature, ios::out);//  sprintf(file_movie, "movie_rho%f_T%f.out", density, reduced_temperature);

  //file velocity
  sprintf(file_velocity, "velocity_rho%f_temperature%f.out", density, reduced_temperature);
  ofstream output_velocity(file_velocity, ios::out);

    //file velocity
    sprintf(file_velocity_all, "velocity_all_rho%f_temperature%f.out", density, reduced_temperature);
    ofstream output_velocity_all(file_velocity_all, ios::out);


  output_energy << 0 << "  " << totalke/particle.size() << "  " << totalpe/particle.size() << "  " << (totalke+totalpe)/particle.size() << endl;

  output_temperature << 0 << "  " << (2.0/3.0)*totalke/particle.size() << endl;

  // print energies
  cout << "initial kinetic energy per particle (in reduced units): " << totalke/particle.size() << endl;
  cout << "initial potential energy per particle (in reduced units): " << totalpe/particle.size() << endl;
  cout << "initial total system energy per particle (in reduced units): " << (totalke+totalpe)/particle.size() << endl;
  cout << "\n";

  double totaltime = 20;
  int steps = 40000;		// number of time discretizations (slices)
  double delta_t = totaltime/steps;	// choose steps carefully
  int movie_step = 100;  // movie will be made every movie_step
  int energycalc_step = 100; // energies will be computed every energycalc_step

  cout << "Fluid simulation time (in seconds): " << totaltime*unittime << " and in reduced units: " << totaltime << endl;

  // Parameters and resources for setting up an NVT ensemble
  double damp = 1;
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::uniform_real_distribution<> get_random_number(-0.5,0.5);

  // code to setup computing averages from simulations
  int hit_eqm = 1000; // this is your choice of where you think the system hit equilibrium
  double average_pe = 0.0;
  double average_ke = 0.0;
  int samples = 0;

  // Molecular Dynamics
  cout << "progress..." << endl;



  double temp_ramp_level = 1; //Change this variable depending on the ramp of temperature change
  int n = 1;
  int k = 1;
  for (int num=1; num <= steps; num++)
  {

//The Ramping code can be read as followed: the variable choose determines how fast you want to ramp the temps.  Each ramp changes temperature after a certain amount
//of steps.  Then, it checks if the number of steps currently is an integer multiple of the number of steps it takes for the temp to change an n number of times
//if it is, it then increases/decreases the temperature an amount that will eventually result in the final temp being reached before the very last step.

      //Rapid Temperature Change - temp changes once
      if(temp_ramp_level==1) {
          if (num == n * 20000) {
              reduced_temperature = reduced_temperature + temperature_change;
              n++;
          }
      }

      //Median Fast Temperature Change - temp changes every 4000 steps, 5 minus 1 times, 4
      if(temp_ramp_level == 4){
          if (num == n*(8000) && num < steps) {
              reduced_temperature = reduced_temperature + (1/temp_ramp_level) * temperature_change;
              n++;
          }
      }

      //Medium Temperature Change - temp changes every 2000 steps -> 10 minus 1 times, 9
      if(temp_ramp_level == 9){
          if (num == n * 4000) {
              reduced_temperature = reduced_temperature + (1 / temp_ramp_level) * temperature_change;
              n++;
          }
      }

      //Median Slow Temperature Change - temp changes every 200 steps -> 100 minus 1 times, 99
      if(temp_ramp_level == 99) {
          if (num == n * 400) {
              reduced_temperature = reduced_temperature + (1 / temp_ramp_level) * temperature_change;
              n++;
          }
      }

      //Slow Temperature Change - temp changes every 10 steps -> 2000 minus 1 times, 1999
      if(temp_ramp_level == 1999) {
          if (num == n * 20) {
              reduced_temperature = reduced_temperature + (1 / temp_ramp_level) * temperature_change;
              n++;
          }
      }

      //Compute the Average velocity - ask how to obtain the average velocities and then should be able to compute magnitude of them from there



      //COMPUTE fdrag
      for(unsigned int i=0; i < particle.size(); i++) {
          particle[i].fdrag = particle[i].velocity * (-particle[i].m * damp);
      }

      // velocity-Verlet
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_velocity(delta_t/2);  // update velocity half timestep
      }

      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_position(delta_t, boxL);  // update position full timestep
      }

      compute_forces(particle, reduced_ljenergy, distance_cutoff, boxL);  // expensive step

      //COMPUTE frandom
      for(unsigned int i=0; i < particle.size(); i++) {
          particle[i].frandom.x = sqrt(24*particle[i].m*reduced_temperature*damp/delta_t)*get_random_number(generator);//particle[i].velocity * (-particle[i].m * damp);
          particle[i].frandom.y = sqrt(24*particle[i].m*reduced_temperature*damp/delta_t)*get_random_number(generator);
          particle[i].frandom.z = sqrt(24*particle[i].m*reduced_temperature*damp/delta_t)*get_random_number(generator);
      }

      //ADD fdrag and frandom to the force computed via compute_forces
      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].force = particle[i].force + particle[i].fdrag + particle[i].frandom;    // you need to add the fdrag and frandom forces on the RHS here
      }

      for (unsigned int i = 0; i < particle.size(); i++) {
          particle[i].update_velocity(delta_t/2);  // update velocity half timestep
      }

      // calling a movie function to get a movie of the simulation every movie_step
      if (num%movie_step == 0)
          make_movie(num,particle,list_propagation,boxL);

      // calculating energies every energycalc_step
      if (num%energycalc_step == 0) {
          totalke = 0.0;
          //calculating avg velo
          VECTOR3D temp_avg_velo = VECTOR3D(0, 0, 0);
          VECTOR3D temp_inst_velo = VECTOR3D(0, 0, 0);

          for(unsigned int  i = 0; i < particle.size(); i++){
              temp_avg_velo = temp_avg_velo + particle[i].velocity;
              temp_inst_velo.x = abs(temp_inst_velo.x) + abs(particle[i].velocity.x);
              temp_inst_velo.y = abs(temp_inst_velo.y) + abs(particle[i].velocity.y);
              temp_inst_velo.z = abs(temp_inst_velo.z) + abs(particle[i].velocity.z);


              //outputting all particles velocity
              //cout << "BALLS BALLS BALLS BALLS" << endl;
              output_velocity_all << num << "  " << particle[i].velocity.x << "  " << particle[i].velocity.y << "  " << particle[i].velocity.z << endl;



          }
          k = k+1;

          VECTOR3D avg_velo = temp_avg_velo*(1/particle.size()); //dividing over number of particles
          cout << "average velocity: " <<  avg_velo.x  <<  " "  <<  avg_velo.y <<  " " <<  avg_velo.z  << endl;

          VECTOR3D inst_velo = temp_inst_velo;
          cout << "instantaneous velocity: " <<  inst_velo.x  <<  " "  <<  inst_velo.y <<  " " <<  inst_velo.z  << endl;

          for (unsigned int i = 0; i < particle.size(); i++) {
              particle[i].kinetic_energy();
              totalke += particle[i].ke;
          }
          totalpe = compute_potential_energies(particle, reduced_ljenergy, distance_cutoff, boxL);

          // outputting the energy (per particle) to make sure simulation can be trusted
          output_energy << num << "  " << totalke/particle.size() << "  " << totalpe/particle.size() << "  " << (totalke+totalpe)/particle.size() << endl;

          // outputting the temperature
          output_temperature << num << "  " << (2.0/3.0)*totalke/particle.size() << endl;


          //outputting_velocity
          output_velocity << num << "  " << inst_velo.x << "  " << inst_velo.y << "  " << inst_velo.z << endl;


          // if equilibrium (steady-state) is reached, we measure the average equilibrium pe and ke by utilizing the generated samples (totalpe, total ke).
          // other equilibrium properties can be measured similarly
          if (num > hit_eqm)
          {
              average_pe = average_pe + totalpe;
              average_ke = average_ke + totalke;
              samples++;
          }
      }

      // monitoring progress of the simulation
      double progress = ((num)/(double)steps);
      ProgressBar(progress);
  }

  cout << endl;
  cout << "equilibrium properties: ..." << endl;
  cout << "average equilibrium pe per particle is " << (average_pe/samples/particle.size()) << endl;
  cout << "average equilibrium ke per particle is " << (average_ke/samples/particle.size()) << endl;

  double equilibrium_temperature = (2.0/3.0)*(average_ke/samples)/(particle.size());
  cout << "average equilibrium temperature is " << equilibrium_temperature << endl;

  return 0;
} 
// End of main