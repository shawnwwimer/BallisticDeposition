#include "SimulationParameters.h"

SimulationParametersFull::~SimulationParametersFull()
{
	clearLayers();
}

void SimulationParametersFull::clearLayers()
{
	for (uint8_t i = 0; i < number_of_layers; i++)
	{
		turns -= parameters[number_of_layers - i - 1]->turns;
		deposited -= parameters[number_of_layers - i - 1]->repetitions;
		delete(parameters[number_of_layers - i - 1]);
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
				serialization += "1";
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
		serialization += "\"Seed\": " + std::to_string(p->seed) + ", ";
		serialization += "\"Time taken\": " + std::to_string(p->time_taken) + ", ";
		serialization += "\"Time finished\": " + std::to_string(p->time_finished) + "}}";
		if (i < number_of_layers - 1) {
			serialization += ", ";
		}
	}
	serialization += "}";
}