//Using standard c++ functions let's create the dictionary with all required variables  
std::ofstream file("myOutputFile", std::ios::app);

//pthread_mutex_t FileMutex = PTHREAD_MUTEX_INITIALIZER;
/*
//Head of the dict
file << "FoamFile" << std::endl << "{"<< std::endl;
file << '\t' << "version 2.0;" << std::endl;
file << '\t' << "format ascii;" << std::endl;
file << '\t' << "class dictionary;" << std::endl;
file << '\t' << "location \"\";" << std::endl; //change location if needed
file << '\t' << "object myChangingDict;" << std::endl; // use your name
file << "}" << std::endl;
*/
// if particles are inside the specified region
if(ii > 1)
{
//pthread_mutex_lock(&FileMutex);

// Number of particles
file << "Time = " << runTime.timeName() << ";" << std::endl;
file << "N" << '\t' << activeList.size() << ";" << std::endl;

//Position
printVectorList(partPos, file, "position", activeList);
//Diameter
printScalarList(partDiam, file, "d", activeList);
//Volume
printScalarList(partVol, file, "V", activeList);
//Velocity
printVectorList(partVel, file, "U", activeList);

//pthread_mutex_unlock(&FileMutex);
}
//else file << "N" << '\t' << 0 << ";" << std::endl;

file.close();
