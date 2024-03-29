#include "SimulationParameters.h"

SimulationParametersFull::~SimulationParametersFull()
{
	clearLayers();
}

void SimulationParametersFull::clearLayers()
{
	int num_layers = number_of_layers;
	for (uint8_t i = 0; i < num_layers; i++)
	{
		turns -= parameters[i]->turns;
		deposited -= parameters[i]->repetitions;
		delete(parameters[i]);
		number_of_layers -= 1;
	}
}

void SimulationParametersFull::addLayer(SimulationParametersIndividual * layer)
{
	parameters[number_of_layers] = layer;
	number_of_layers += 1;
	turns += layer->turns;
	deposited += layer->repetitions;
}

void SimulationParametersFull::serialize() {
	// Stack parameters
	serialization = "{\"Layers\": " + std::to_string(number_of_layers) + ", \"Total points\": " + std::to_string(deposited) + ", ";

	// Go through each parameter
	for (int i = 0; i < number_of_layers; i++) {
		SimulationParametersIndividual* p = parameters[i];
		std::string theta_str = std::to_string(p->theta);
		theta_str.erase(theta_str.find_last_not_of('0') + 1, std::string::npos);
		theta_str.erase(theta_str.find_last_not_of('.') + 1, std::string::npos);
		std::string phi_str = std::to_string(p->phi);
		phi_str.erase(phi_str.find_last_not_of('0') + 1, std::string::npos);
		phi_str.erase(phi_str.find_last_not_of('.') + 1, std::string::npos);
		std::string turns_str = std::to_string(p->turns);
		turns_str.erase(turns_str.find_last_not_of('0') + 1, std::string::npos);
		turns_str.erase(turns_str.find_last_not_of('.') + 1, std::string::npos);
		std::string phi_deg_str = std::to_string(p->phi_deg);
		phi_deg_str.erase(phi_deg_str.find_last_not_of('0') + 1, std::string::npos);
		phi_deg_str.erase(phi_deg_str.find_last_not_of('.') + 1, std::string::npos);
		std::string theta_end_str = std::to_string(p->theta_end);
		theta_end_str.erase(theta_end_str.find_last_not_of('0') + 1, std::string::npos);
		theta_end_str.erase(theta_end_str.find_last_not_of('.') + 1, std::string::npos);
		serialization += "\"" + std::to_string(i + 1) + "\": {";
		serialization += "\"Deposited\": " + std::to_string(p->deposited) + ", ";
		serialization += "\"Parameters\": {";
		//serialization += "\"System\": " + p->system + ", ";
		serialization += "\"L\": " + std::to_string(p->length) + ", ";
		serialization += "\"W\": " + std::to_string(p->width) + ", ";
		serialization += "\"theta\": " + theta_str + ", ";
		serialization += "\"phi\": " + phi_str + ", ";
		serialization += "\"H\": " + std::to_string(p->height) + ", ";
		serialization += "\"D\": " + std::to_string(p->diffusion_steps) + ", ";
		serialization += "\"turns\": " + turns_str + ", ";
		serialization += "\"phi sweeps\": " + std::to_string(p->phi_num) + ", ";
		serialization += "\"phi sweep degrees\": " + phi_deg_str + ", ";
		serialization += "\"theta end\": " + theta_end_str + ", ";
		serialization += "\"stepper resolution\": " + std::to_string(p->stepper_resolution) + ", ";
		serialization += "\"species\": [";
		for (int s = 0; s < p->species->size(); s++) {
			std::string str = std::to_string((*(p->species))[s]);
			str.erase(str.find_last_not_of('0') + 1, std::string::npos);
			str.erase(str.find_last_not_of('.') + 1, std::string::npos);
			serialization += str;
			if (s < p->species->size() - 1) {
				serialization += ", ";
			}
			else {
				serialization += "], ";
			}
		}
		serialization += "\"weights\": [";
		for (int w1 = 0; w1 < p->weights->size(); w1++) {
			serialization += "[";
			for (int w2 = 0; w2 < p->weights->size(); w2++) {
				serialization += std::to_string((*(p->weights))[w1][w2]);
				if (w2 < p->weights->size() - 1) {
					serialization += ", ";
				}
				else {
					if (w1 < p->weights->size() - 1) {
						serialization += "], ";
					}
					else {
						serialization += "]";
					}
					
				}
			}
		}
		serialization += "], ";
		serialization += "\"repetition\": " + std::to_string(p->repetitions) + ", ";
		serialization += "\"spread\": [";
		for (int s = 0; s < p->spread->size(); s++) {
			std::string str = std::to_string((*(p->spread))[s]);
			str.erase(str.find_last_not_of('0') + 1, std::string::npos);
			str.erase(str.find_last_not_of('.') + 1, std::string::npos);
			serialization += str;
			if (s < p->spread->size() - 1) {
				serialization += ", ";
			}
			else {
				serialization += "], ";
			}
		}
		if (p->acceleration == 0) {
			serialization += "\"Acceleration method\": NONE, ";
		}
		else if (p->acceleration == 1) {
			serialization += "\"Acceleration method\": ACC, ";
		}
		else if (p->acceleration == 2) {
			serialization += "\"Acceleration method\": DEC, ";
		}
		else if (p->acceleration == 3) {
			serialization += "\"Acceleration method\": BICONE, ";
		}
		else if (p->acceleration == 4) {
			serialization += "\"Acceleration method\": HOURGLASS, ";
		}
		else {
			serialization += "\"Acceleration method\": undefined, ";
		}
		serialization += "\"Acceleration method code\": " + std::to_string(p->acceleration) + ", ";

		if (p->collision_method == 0) {
			serialization += "\"Collision method\": NNO, ";
		}
		else if(p->collision_method == 1) {
			serialization += "\"Collision method\": NN1, ";
		}
		else if(p->collision_method == 2) {
			serialization += "\"Collision method\": NN2, ";
		}
		else if(p->collision_method == 3) {
			serialization += "\"Collision method\": NN3, ";
		}
		else {
			serialization += "\"Collision method\": undefined, ";
		}
		serialization += "\"Acceleration method code\": " + std::to_string(p->collision_method) + ", ";
		
		serialization += "\"Seed\": " + std::to_string(p->seed) + ", ";
		serialization += "\"Time taken\": " + std::to_string(p->time_taken) + ", ";
		serialization += "\"Time finished\": " + std::to_string(p->time_finished) + "}}";
		if (i < number_of_layers - 1) {
			serialization += ", ";
		}
	}
	serialization += "}";
}


ContinuousSimulationParametersFull::~ContinuousSimulationParametersFull()
{
	clearLayers();
}

void ContinuousSimulationParametersFull::clearLayers()
{
	int num_layers = number_of_layers;
	for (uint8_t i = 0; i < num_layers; i++)
	{
		deposited -= parameters[i]->repetitions;
		delete(parameters[i]);
		number_of_layers -= 1;
	}
}

void ContinuousSimulationParametersFull::addLayer(ContinuousSimulationParameters* layer)
{
	parameters[number_of_layers] = layer;
	number_of_layers += 1;
	deposited += layer->repetitions;
}

void ContinuousSimulationParametersFull::serialize() {
	// Stack parameters
	serialization = "{\"Layers\": " + std::to_string(number_of_layers) + ", \"Total points\": " + std::to_string(deposited) + ", ";

	// Go through each parameter
	for (int i = 0; i < number_of_layers; i++) {
		ContinuousSimulationParameters* p = parameters[i];
		std::string theta_str = std::to_string(p->theta);
		theta_str.erase(theta_str.find_last_not_of('0') + 1, std::string::npos);
		theta_str.erase(theta_str.find_last_not_of('.') + 1, std::string::npos);
		std::string diff_str = std::to_string(p->diffusion_length);
		diff_str.erase(diff_str.find_last_not_of('0') + 1, std::string::npos);
		diff_str.erase(diff_str.find_last_not_of('.') + 1, std::string::npos);
		serialization += "\"" + std::to_string(i + 1) + "\": {";
		serialization += "\"Repetitions\": " + std::to_string(p->repetitions) + ", ";
		serialization += "\"Parameters\": {";
		//serialization += "\"System\": " + p->system + ", ";
		serialization += "\"L\": " + std::to_string(p->length) + ", ";
		serialization += "\"W\": " + std::to_string(p->width) + ", ";
		serialization += "\"theta\": " + theta_str + ", ";
		serialization += "\"H\": " + std::to_string(p->height) + ", ";
		serialization += "\"D\": " + diff_str + ", ";
		serialization += "\"species\": [";
		for (int s = 0; s < p->species->size(); s++) {
			std::string str = std::to_string((*(p->species))[s]);
			str.erase(str.find_last_not_of('0') + 1, std::string::npos);
			str.erase(str.find_last_not_of('.') + 1, std::string::npos);
			serialization += str;
			if (s < p->species->size() - 1) {
				serialization += ", ";
			}
			else {
				serialization += "], ";
			}
		}
		serialization += "\"radii\": [";
		for (int r = 0; r < p->radii->size(); r++) {
			std::string str = std::to_string((*(p->radii))[r]);
			str.erase(str.find_last_not_of('0') + 1, std::string::npos);
			str.erase(str.find_last_not_of('.') + 1, std::string::npos);
			serialization += str;
			if (r < p->radii->size() - 1) {
				serialization += ", ";
			}
			else {
				serialization += "], ";
			}
		}
		serialization += "\"bin size\": " + std::to_string(p->bin_size) + ", ";
		std::string scale_str = std::to_string(p->length_scale);
		scale_str.erase(scale_str.find_last_not_of('0') + 1, std::string::npos);
		scale_str.erase(scale_str.find_last_not_of('.') + 1, std::string::npos);
		serialization += "\"length scale\": " + scale_str + ", ";
		std::string cube_str = std::to_string(p->cube_size);
		cube_str.erase(cube_str.find_last_not_of('0') + 1, std::string::npos);
		cube_str.erase(cube_str.find_last_not_of('.') + 1, std::string::npos);
		serialization += "\"cube size\": " + cube_str + ", ";
		serialization += "\"diffusion method\": " + std::to_string(p->diffusion_method) + ", ";
		serialization += "\"repetition\": " + std::to_string(p->repetitions) + ", ";
		serialization += "\"Seed\": " + std::to_string(p->seed) + ", ";
		serialization += "\"Time taken\": " + std::to_string(p->time_taken) + ", ";
		serialization += "\"Time finished\": " + std::to_string(p->time_finished) + "}}";
		if (i < number_of_layers - 1) {
			serialization += ", ";
		}
	}
	serialization += "}";
}