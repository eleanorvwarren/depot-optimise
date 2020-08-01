// main.cpp
// a program that reads in the csv data file GBPlaces and computes optimised locations for two delivery hubs
// by minimising the total distance from each hub to the places it serves, weighted by the population of each place
// where the places each hub serves are selected by choosing the hub closest to the place

// Eleanor Warren 12/2017

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#define PI 3.14159

using namespace std;

// define a structure for GBPlaces data
struct GBPlace {
    string place, type;
    double population, latitude, longitude;
};

// define a structure for holding information for hub 1 and hub 2 simultaneously
struct twoHubs {
    double hub1, hub2;
};


// function definitions

// random number generator
double random_number ( double upper, double lower, int n ) {
    // function to return a random number between two limits, lower and upper
    // n is the number of bits into which to split the range
    double r;
    r = lower + (rand() % (n + 1) * (1./n) * (upper-lower));
    return r;
}

// function to calculate the great circle distance between two points
double haversine (double lat1, double lat2, double long1, double long2) {
    // convert latitudes and longitudes to radians
    lat1 = lat1*(PI/180);
    lat2 = lat2*(PI/180);
    long1 = long1*(PI/180);
    long2 = long2*(PI/180);

    // apply the haversine formula
    double R = 3958.75;
    double dLat = lat2 - lat1;
    double dLong = long2 - long1;
    double a = pow(sin(dLat/2),2) + cos(lat1)*cos(lat2)*sin(dLong/2)*sin(dLong/2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    double distance = R * c;

    return distance;
}

// fitness function
// function that compares the distance from each hub to one place, and adds the smaller distance to the corresponding hubs total distance
twoHubs totalDistanceDistributor (vector <GBPlace> GBPlaces_data, double trial_lat_1, double trial_long_1, double trial_lat_2, double trial_long_2) {
    // define variables
    double haversine_1, haversine_2;
    double total_distance_1, total_distance_2;
    // loop through all places
    for (int i = 0; i < GBPlaces_data.size(); i++) {
        // calculate the distance from the trial lat,long of each hub to each place, and weight by the sqrt of its population
        haversine_1 = haversine(trial_lat_1, GBPlaces_data[i].latitude, trial_long_1, GBPlaces_data[i].longitude)/sqrt(GBPlaces_data[i].population);
        haversine_2 = haversine(trial_lat_2, GBPlaces_data[i].latitude, trial_long_2, GBPlaces_data[i].longitude)/sqrt(GBPlaces_data[i].population);
        // and add it onto the sum variable for the hub with the lowest distance
        if (haversine_1 < haversine_2) {
            total_distance_1 += haversine_1;
        } else {
            total_distance_2 += haversine_2;
        }

    }
    // output the results as a twoHubs structure of the total distances
    twoHubs total_distances;
    total_distances.hub1 = total_distance_1;
    total_distances.hub2 = total_distance_2;

    return total_distances;
}

// function to calculate the magnitude of values of structure twoHubs
double magnitude (twoHubs valuesSet) {
    double magnitude = valuesSet.hub1 + valuesSet.hub2;
    return magnitude;
}

// main program
int main() {
    // declare variables
    int d_trial_lat_1, d_trial_long_1, d_trial_lat_2, d_trial_long_2; // holds the most recent change in latitude and longitude positions of the two hubs
    double trial_lat_1, trial_long_1, trial_lat_2, trial_long_2; // holds the current best values of latitude and longitude for each hub
    double step = 0.01; // size of step to move in latitude and longitude
    twoHubs values, oldValues, newValues, newValues1, newValues2, newValues12;
    double global_min = 1e10; // large starting value of global min to be decreased later on
    double global_min_lat_1, global_min_long_1, global_min_lat_2, global_min_long_2; // coordinates for each of the hubs in the optimised positions
    string line;
    vector <GBPlace> GBPlaces_data; // vector of structures of type GBPlace
    GBPlace place_data; // structure to store the data for a place temporarily

    srand(time(NULL)); // seeds random number generator

    // read in GBplaces.csv
    ifstream theFile ("../GBplaces.csv");

    // declare variables for data in theFile
    string place, type, population, latitude, longitude;

    getline(theFile, line);

    while( getline(theFile, place,',') ) {
        // separate the file into its data components
        getline(theFile, type, ',');
        getline(theFile, population, ',');
        getline(theFile, latitude, ',');
        getline(theFile, longitude);

        // add the data to the temporary GBPlace structure
        place_data.place = place;
        place_data.type = type;
        place_data.population = atof(population.c_str());
        place_data.latitude = atof(latitude.c_str());
        place_data.longitude = atof(longitude.c_str());

        // add the data to the vector of structures
        GBPlaces_data.push_back(place_data);

    }

    theFile.close(); // close the file


    for (int min_counter = 0; min_counter < 100; min_counter ++) { // run the program multiple times from random starting points to get the local minima

        // pick a random starting point for latitude and longitude of each of the two hubs
        trial_lat_1 = random_number(50.3, 57.2, 100);
        trial_lat_2 = random_number(50.3, 57.2, 100);
        trial_long_1 = random_number(-4.3, -1.3, 100);
        trial_long_2 = random_number(-4.3, -1.3, 100);

        // calculate values of the function totalDistanceDistributor for the two hub locations
        values = totalDistanceDistributor (GBPlaces_data, trial_lat_1, trial_long_1, trial_lat_2, trial_long_2);

        // print the current values
        cout << "Trial " << min_counter << ": \n" << values.hub1 << " & " << values.hub2 << "\n";
        // find a local minimum value
        do {
            // look around the current positions and see if we can move to decrease the values
            oldValues = values; // redefine the previous values as oldValues


            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    if (i == 0 && j == 0) {
                        // position has not moved, do nothing
                    } else {
                        // calculate the magnitudes of the value if hub 1 and/or hub2 move by i and j respective to their positions
                        // this allows the hubs to move independently of one another
                        newValues1 = totalDistanceDistributor(GBPlaces_data, trial_lat_1 + step*i, trial_long_1 + step*j, trial_lat_2, trial_long_2);
                        newValues2 = totalDistanceDistributor(GBPlaces_data, trial_lat_1, trial_long_1, trial_lat_2 + step*i, trial_long_2 + step*j);
                        newValues12 = totalDistanceDistributor(GBPlaces_data, trial_lat_1 + step*i, trial_long_1 + step*j, trial_lat_2 + step*i, trial_long_2 + step*j);

                        // select the option out of the three that gives the lowest sum of total distances
                        // define a variable to keep track of the choice made
                        int c;

                        if (magnitude(newValues1) < magnitude(newValues2) && magnitude(newValues1)<magnitude(newValues12)){
                            newValues = newValues1;
                            c=1;
                        } else if (magnitude(newValues2) < magnitude(newValues1) && magnitude(newValues2) < magnitude(newValues12)) {
                            newValues = newValues2;
                            c=2;
                        } else {
                            newValues = newValues12;
                            c=12;
                        }

                        // if this new value is less than the previous value, rewrite it as the value and save it's relative position from the previous point
                        if (magnitude(newValues) <= magnitude(values)) {
                            values = newValues;
                            // use the variable c to add the i and j to the correct hub(s)
                            if (c==1) {
                                d_trial_lat_1 = i;
                                d_trial_long_1 = j;
                                d_trial_lat_2 = 0;
                                d_trial_long_2 = 0;
                            } else if (c==2){
                                d_trial_lat_1 = 0;
                                d_trial_long_1 = 0;
                                d_trial_lat_2 = i;
                                d_trial_long_2 = j;
                            } else if (c==12) {
                                d_trial_lat_1 = i;
                                d_trial_long_1 = j;
                                d_trial_lat_2 = i;
                                d_trial_long_2 = j;
                            }

                        }
                    }
                }
            }

            // update the hub positions and the values by adding on the chosen changes in lat and long
            trial_lat_1 += step * d_trial_lat_1;
            trial_long_1 += step * d_trial_long_1;
            trial_lat_2 += step * d_trial_lat_2;
            trial_long_2 += step * d_trial_long_2;

            // print values as they are computed
            cout << values.hub1 << " & " << values.hub2 << "\n";

        } while (magnitude(values) < magnitude(oldValues)); // run through this loop while the value is still decreasing

        // print the local minima found for the initial random starting positions
        cout << "Local minima " << min_counter << " at (" << trial_lat_1 << "," << trial_long_1 << ") and (" << trial_lat_2 << "," << trial_long_2 << ") give total value " << magnitude(values) << "\n\n";

        // if the value for the local minima found is less than the value for global minima, make it the new global minima value
        if (magnitude(values) < global_min){
            global_min = magnitude(values);
            // make the corresponding coordinates the global_min coordinates
            global_min_lat_1 = trial_lat_1;
            global_min_long_1 = trial_long_1;
            global_min_lat_2 = trial_lat_2;
            global_min_long_2 = trial_long_2;
        }
    }

    // output result for the global minima, ie the optimised locations for the two hubs
    cout << "The optimal locations of two delivery hubs are at coordinates (" << global_min_lat_1 << "," << global_min_long_1 << ") and (" << global_min_lat_2 << "," << global_min_long_2 << ").";

    return 0;
}
