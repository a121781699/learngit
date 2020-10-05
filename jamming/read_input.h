#ifndef READ_INPUT_H
#define READ_INPUT_H

#define NAME_LENGTH 256

class read_input
{
public:
	read_input() {}
	//box parameter
	int N;
	int ignore_line;
	char initial_configuration[NAME_LENGTH];
	char final_configuration[NAME_LENGTH];
	bool monitor;
	char MonitorFile[NAME_LENGTH];
	int read(int argc, char* argv[]);

	//FIRE parameter
	double phi;
	int Nmin;
	int read_FIRE(int argc, char* argv[]);

	~read_input();
};
#endif
