In my program I have chosen to optimise the locations of two delivery hubs by minimising the total distance from a hub to its assigned places. I used the assumption that the delivery vans drive from the hubs outwards to a place and back before another delivery. I have assigned each place in the GBplaces.csv file to one of the hubs, based on which great circle distance weighted by the population of that place is smaller, so that hubs are also closer to the places where more deliveries are likely to be necessary. The locations of each hub can move independently of one another, in order for both to be optimised. The results for all local minima computed are clearly printed, as well as the final solutions at the end. 