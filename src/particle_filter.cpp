/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	Particle single_particle;

	//set number of particles
	num_particles = 100;

	// This line creates a normal (Gaussian) distribution for x.
	normal_distribution<double> dist_x(x, std[0]);

	//Create normal distributions for y and psi.

	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);

	for (int i = 0; i < num_particles; ++i) {

		// Sample from these normal distrubtions like this:
		// sample_x = dist_x(gen);
		// where "gen" is the random engine initialized earlier .
		single_particle.id = i+1;
		single_particle.x = dist_x(gen);
		single_particle.y = dist_y(gen);
		single_particle.theta = dist_psi(gen);
		single_particle.weight = 1.0;
		particles.push_back(single_particle);
	}

	is_initialized = true;

}



void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	//default_random_engine gen;

	double x_pred_mean;
	double y_pred_mean;
	double theta_mean;

	for (int i = 0; i < num_particles; i++) {


		//calculate predicted x, y, theta

		if (fabs(yaw_rate) < 0.00001) {
			x_pred_mean = particles[i].x + velocity * delta_t * cos(particles[i].theta);
			y_pred_mean = particles[i].y + velocity * delta_t * sin(particles[i].theta);
		}
		else {
			x_pred_mean = particles[i].x + (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			y_pred_mean = particles[i].y + (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
			theta_mean = particles[i].theta + yaw_rate * delta_t;
		}

		//define generator with predicted x, y, theta as mean and with std_pos noise
		normal_distribution<double> dist_x(x_pred_mean, std_pos[0]);
		normal_distribution<double> dist_y(y_pred_mean, std_pos[1]);
		normal_distribution<double> dist_psi(theta_mean, std_pos[2]);

		//predict particle with noise
		particles[i].x = dist_x(gen);
		particles[i].y = dist_y(gen);
		particles[i].theta = dist_psi(gen);

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	double distance;

	for (int k = 0; k < observations.size(); k++) {

		//calculate initial distance according to the first landmark item and set initial observation ID to the first landmark id
		distance = dist(observations[k].x, observations[k].y, predicted[0].x, predicted[0].y);
		observations[k].id = predicted[0].id;

		//looping through landmark items and updating current observation id if distance is smaller
		for (int i = 1; i < predicted.size(); i++) {
			if (dist(observations[k].x, observations[k].y, predicted[i].x, predicted[i].y) < distance) {
				distance = dist(observations[k].x, observations[k].y, predicted[i].x, predicted[i].y);
				observations[k].id = predicted[i].id;
			}

		}

	}

}



void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html



	vector<LandmarkObs> landmarks_in_range;
	vector<LandmarkObs> trans_observations;

	LandmarkObs single_observation;

	for (int i = 0; i < num_particles; i++) {

		landmarks_in_range.clear();
		trans_observations.clear();

		//keep landmarks within sensor range
		for (int n = 0; n < map_landmarks.landmark_list.size(); n++) {
			if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[n].x_f, map_landmarks.landmark_list[n].y_f) <= sensor_range) {
				landmarks_in_range.push_back(LandmarkObs{map_landmarks.landmark_list[n].id_i, map_landmarks.landmark_list[n].x_f, map_landmarks.landmark_list[n].y_f});
			}
		}

		//transform observations from car to map coordinate
		for (int j = 0; j < observations.size(); j++) {

			double t_x = particles[i].x + observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta);
			double t_y = particles[i].y + observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta);

			trans_observations.push_back(LandmarkObs{-1, t_x, t_y});
		}

		//get landmark id for each transformed observation
		dataAssociation(landmarks_in_range, trans_observations);

		//calculate weight for each particle based on assigned transformed observations and landmark coordinates
		for (int m = 0; m < trans_observations.size(); m++) {

			double lm_x;
			double lm_y;

			for (int n = 0; n < landmarks_in_range.size(); n++) {
				if (landmarks_in_range[n].id == trans_observations[m].id) {
					lm_x = landmarks_in_range[n].x;
					lm_y = landmarks_in_range[n].y;
				}
			}

			double dev_x = pow((trans_observations[m].x - lm_x),2);
			double dev_y = pow((trans_observations[m].y - lm_y),2);
			particles[i].weight *= 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]) * exp (-(dev_x / (2 * pow(std_landmark[0],2)) + (dev_y / (2 * pow(std_landmark[1],2)))));
		}

	}

}



void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	//default_random_engine generator;
	//discrete_distribution<int> disc_dist(0,num_particles-1);
	uniform_int_distribution<int> disc_dist(0, num_particles-1);

	int index = disc_dist(gen);

	vector<Particle> resamp_particles;
	Particle single_particle;

	weights.clear();

	double sum = 0.0;
	double beta = 0.0;

	//calculate normalized weight
	for (int i = 0; i < num_particles; i++) {
		sum += particles[i].weight;
	}

	for (int i = 0; i < num_particles; i++) {
		weights.push_back(particles[i].weight/sum);
	}

	double max_weight = *max_element(weights.begin(), weights.end());

	//uniform generator;
	uniform_real_distribution<double> unidistribution(0.0, max_weight * 2);

	//resampling
	for (int i = 0; i < num_particles; i++) {
		beta = beta + unidistribution(gen);

		while (weights[index] < beta) {
			beta = beta - weights[index];
			if (index == num_particles) {
				index = 0;
			} else {
				index++;
			}
		}

		single_particle.x = particles[index].x;
		single_particle.y = particles[index].y;
		single_particle.theta = particles[index].theta;
		single_particle.weight = 1.0;

		resamp_particles.push_back(single_particle);

	}

	//save new particles
	particles = resamp_particles;


}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
