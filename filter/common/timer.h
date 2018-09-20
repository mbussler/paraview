/* timer writeout */
#include <stdlib.h> /* getenv */
#include <fstream>
#include <string>

void write_timer( std::string classname, std::string functionname, int time )
{
    // performance measurement
    char* logFilename = NULL;
    logFilename = getenv("PV_LOGFILE");
    if( logFilename != NULL )
    {
        std::string filenname(logFilename);
        filenname += "_" + classname + "_" + functionname + ".txt";
        std::ofstream logFile(filenname.c_str(), std::ios_base::app);

        logFile << time << "\n";
    }

    std::cout << "[" << classname << "]" << "[" << functionname << "] " << time << " ms\n";
}
