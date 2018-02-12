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
#define EPS 0.00001
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    if (is_initialized == TRUE)
    {return;}

    num_particles = 100;
    double std_x = std[0];
    double std_y = std[1];
    double std_theta= std[2];

    normal_distribution<double> x_dist{x,std_x};
    normal_distribution<double> y_dist{y,std_y};
    normal_distribution<double> theta_dist{theta,std_theta};

    for(i=0;i<num_particles;i++){
        Particle particle;
        particle.id = i;
        particle.x=x_dist(gen);
        particle.y=y_dist(gen);
        particle.theta=theta_dist(gen);
        particle.weight=1.0;
        particles.push_back(particle);
    }
    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
    double std_x=std_pos[0];
    double std_y=std_pos[1];
    double std_theta=std_pos[2];

    normal_distribution<double> x_dist{0,std_x};
    normal_distribution<double> y_dist{0,std_y};
    normal_distribution<double> theta_dist{0,std_theta};


    for (int i=0;i<num_particles;i++){

        if(fabs(yaw_rate)> EPS)
        {
            particles[i].x+= (velocity/yaw_rate)*(sin(theta+yaw_rate*delta_t)-sin(theta));
            particles[i].y+=(velocity/yaw_rate)*(cos(theta)-cos(theta+yaw_rate*delta_t));
            particles[i].theta+=yaw_rate*delta_t;
        }
        else
        {
            particles[i].x+=velocity*delta_t*cos(particles[i].theta);
            particles[i].y+=velocity*delta_t*sin(particles[i].theta);
        }
        particles[i].x+=x_dist(gen);
        particles[i].y+=y_dist(gen);
        particles[i].theta+=theta_dist(gen);
    }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
    for(unsigned int i=0; i<observations.size();i++)
    {
        double minDistance = -1;
        int IdFound = -1;

        for(unsigned int j=0; j<predicted.size();j++)
        {
            double distance = pow((observations[i].x - predicted[j].x),2) +
                    pow((observations[i].y - predicted[j].y),2);

            if(minDistance != -1)
            {
                if(distance < minDistance){
                    minDistance =distance;
                    IdFound = predicted[j].id;
                }
            }else
            {
                minDistance = distance;
                IdFound = predicted[j].id;
            }
        }
        observations[i].id=IdFound;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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
    double stdRange=std_landmark[0];
    double stdBearing = std_landmark[1];
    for (int i=0; i<num_particles;i++) {

        double x =particles[i].x;
        double y=particles[i].y;
        double theta = particles[i].theta;
        vector<LandmarkObs> relLandmarks;
        for (unsigned int j=0; j<map_landmarks.landmark_list.size();j++){
            float landmarkX=map_landmarks.landmark_list[j].x_f;
            float landmarkY=map_landmarks.landmark_list[j].y_f;
            int landmarkId=map_landmarks.landmark_list[j].id_i;
            double dX = x-landmarkX;
            double dY = y-landmarkY;
            if(dX*dX+dY*dY <= sensor_range*sensor_range){
                relLandmarks.push_back(LandmarkObs{landmarkId,landmarkX,landmarkY});

            }
        }
        vector<LandmarkObs> transformedObs;
        for(unsigned int j=0;j<observations.size();j++){
            double x_t = con(theta)*observations[j].x -sin(theta)*observations[j].y +x;
            double y_t = sin(theta)*observations[j].x - cos(theta)*observations[j].y + y;
            transformedObs.push_back(LandmarkObs{observations[j].id, x_t,y_t});
        }
        dataAssociation(relLandmarks, transformedObs);
        particles[i].weight =1.0;
        for(unsigned int j=0;j<transformedObs.size();j++){
            double landmarkX, landmarkY;
            unsigned int k=0;
            bool found = false;
            while(!found && k<relLandmarks.size()){
                if(relLandmarks[k].id == transformedObs[j].id ){
                    found = true;
                    landmarkX = relLandmarks[k].x;
                    landmarkY =relLandmarks[k].y;
                }
                k++;
            }

            double dX= transformedObs[j].x - landmarkX;
            double dY= transformedObs[j].y- landmarkY;

            double weight = (1/(2*M_PI*stdRange*stdBearing))* exp( -( dX*dX/(2*stdLandmarkRange*stdLandmarkRange) + (dY*dY/(2*stdLandmarkBearing*stdLandmarkBearing)) ) );
            if(weight == 0){
                particles[i].weight*=EPS;
            }else{
                particles[i].weight*=weight;
            }
        }
    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

    vector<double> weights;
    double maxWeight =-1;
    //find out the maximum weight
    for(unsigned int i=0;i<num_particles;i++){
        weights.push_back(particles[i].weight);
                if(particles[i].weight>maxWeight){
            maxWeight =particles[i].weight;
        }
    }

    uniform_real_distribution<double> distDouble(0.0, maxWeight);
    uniform_int_distribution<int> distInt(0,num_particles -1);

    //resample the weights
    int index=distInt(gen);
    double beta = 0.0;

    vector<Particle> resampledParticles;
    for(unsigned int i =0; i< num_particles; i++){
        beta+=distDouble(gen)*2.0;
        while(beta>weights[index]){
            beta-=weights[index];
            index=(index+1)%num_particles;
        }
        resampledParticles.push_back(particles[index]);

    }
    particles = resampledParticles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
