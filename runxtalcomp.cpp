//===================================================================
// 
//
//
//===================================================================

#include "xtalcomp.h"

#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#define PRINT_DIV \
  printf("|-------------------------------------|-------------------------------------|\n")

std::string debug;

float StringToDouble(const std::string & instring)
{
    double i;
    std::stringstream ss;
    ss << instring;
    ss >> i;
    return i;
}

bool parsePOSCAR(std::string, XcMatrix&, std::vector<XcVector>&, std::vector<unsigned int>&);

void Debug(const char *str, const double d)
{
  char buffer[128];
  snprintf(buffer, 32, "%s %f\n", str, d);
  debug += buffer;
}
void Debug(const std::string &str, const double d) {Debug(str.c_str(), d);}

bool parsePOSCAR(std::string str, XcMatrix &cell,
                 std::vector<XcVector> &pos,
                 std::vector<unsigned int> &types)
{
    std::string line;
  bool cart = false;

  std::ifstream lines (str.c_str());
    if(!lines.is_open()){
        std::cout << "File does not exist!" << std::endl;
        return 0;
    } else{

        // First line is comment
        getline(lines, line);

  // Next line is scale factor
  getline(lines, line);
  float scale;
  if (sscanf(line.c_str(), "%f", &scale) != 1) return false;

  // Next comes the matrix
  float x,y,z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[0][0] = x;
  cell[0][1] = y;
  cell[0][2] = z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[1][0] = x;
  cell[1][1] = y;
  cell[1][2] = z;
  getline(lines, line);
  if (sscanf(line.c_str(), "%f %f %f",
             &x, &y, &z) != 3) return false;
  cell[2][0] = x;
  cell[2][1] = y;
  cell[2][2] = z;

  // Apply scale:
  cell *= scale;

  // Store frac->cart matrix
  XcMatrix toCart = cell.transpose().inverse();

  // List of atom types
  std::vector<int> counts (15); // Allow up to 15 atom types.
  getline(lines, line);
  int tmpint;
  int numTypes = sscanf(line.c_str(), "%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                        &counts[0], &counts[1], &counts[2],
                        &counts[3], &counts[4], &counts[5],
                        &counts[6], &counts[7], &counts[8],
                        &counts[9], &counts[10], &counts[11],
                        &counts[12], &counts[13], &counts[14], &tmpint);

  if (numTypes > 15) return false;

  // Starts with either [Ss]elective dynamics, [KkCc]artesian, or
  // other for fractional coords.
  getline(lines, line);

  // If selective dynamics, get the next line
  if (line.at(0) == 'S' || line.at(0) == 's')
    getline(lines, line);

  // Check if we're using cartesian or fractional coordinates:
  if (line.at(0) == 'K' || line.at(0) == 'k' ||
      line.at(0) == 'C' || line.at(0) == 'c' )
    cart = true;
  else
    cart = false;


  // Coordinates
  // determine number of atoms:
  types.clear();
  int numAtoms = 0;
  for (int i = 0; i < numTypes; ++i) {
    numAtoms += counts[i];
    for (int j = 0; j < counts[i]; ++j) {
      types.push_back(i);
    }
  }

  types.resize(numAtoms);

  Debug("numAtoms:", numAtoms);

  // Grab vectors
  XcVector tmp;
  pos.clear();
  for (int atom_i = 0; atom_i < numAtoms; ++atom_i) {
    getline(lines, line);
    if (sscanf(line.c_str(), "%f %f %f",
               &x, &y, &z) != 3) return false;
    tmp = XcVector(x,y,z);
    debug += "pos line: " + line + "\n";
    Debug("x: ", tmp.x());
    Debug("y: ", tmp.y());
    Debug("z: ", tmp.z());
    if (cart) {
      tmp = toCart * tmp;
      debug += "Converted to cartesian:\n";
      Debug("x: ", tmp.x());
      Debug("y: ", tmp.y());
      Debug("z: ", tmp.z());
    }
    pos.push_back(tmp);
  }

  Debug("pos size: ", pos.size());
  Debug("types size: ", types.size());
    }
  lines.close();
  return true;
}

int main(int argc, char ** argv) {

  XcMatrix cell1 (4.5, 0.0, 0.0,
                  1.2, 2.4, 0.0,
                  2.5, 6.4, 1.1);

  XcMatrix cell2 (2.5, 0.0, 0.0,
                  2.3, 4.2, 0.0,
                  1.5, 2.1, 2.1);

  std::vector<XcVector> pos1 (5);
  std::vector<XcVector> pos2 (5);

  std::vector<unsigned int> types1 (5);
  std::vector<unsigned int> types2 (5);

  double cartTol = 0.05;
  double angleTol = 0.25;
  float transform[16];
  bool validInput = true;
  std::string poscar1;
  std::string poscar2;
  bool foundPOS1 = false;
  bool foundPOS2 = false;

  for (int i=1; i < argc; i++){
      std::string argstr(argv[i]);
        
      if (argstr=="-h" || argstr=="--help"){
          std::cout << std::endl;
          std::cout << "   XtalComp" << std::endl
            << "   Usage:   xtalcomp [options] <POSCAR1> <POSCAR2>" << std::endl << std::endl
            << "   Info:    Input the POSCARs of the two structures you want to compare." << std::endl 
            << "            Enter the tolerances you want (or use the defaults)." << std::endl 
            << "            The results will be printed to the command line." << std::endl
            << std::endl
            << "   Options:" << std::endl
            << "    -h or --help                            Help. This is it!" << std::endl
            << "    -c or --cart <double>                   Determines the cartesian tolerance" << std::endl
            << "                                               Default is 0.05 angstroms." << std::endl
            << "    -a or --ang <double>                    Determines the angle tolerance" << std::endl
            << "                                               Default is 0.25 degrees." << std::endl;
        std::cout << std::endl;
        return 0;
                                                                                                                                     }


    else if (argstr=="-c" || argstr=="--cart"){
        cartTol = StringToDouble(argv[++i]);
        //if (sscanf(argv[++i], "%f", &cartTol) != 1) validInput = false;
        //std::cout << "cart is " << cartTol << "\n";
    }

    else if (argstr=="-a" || argstr=="--ang"){
        angleTol = StringToDouble(argv[++i]);
        //if (sscanf(argv[++i], "%f", &angleTol) != 1) validInput = false;
        //std::cout << "angle is " << angleTol << "\n";
    }

    else{

        if (!foundPOS1){
            poscar1 = argstr;
            foundPOS1 = true;
            if (!parsePOSCAR(poscar1, cell1, pos1, types1)) validInput = false;
        }
        
        else {
            poscar2 = argstr;
            foundPOS2 = true;
            if (!parsePOSCAR(poscar2, cell2, pos2, types2)) validInput = false;
        }
    }
  }

  if (!validInput) {
    printf("XtalComp Results\n");
    printf("Invalid input\n");
    printf("Go back and check your inputs.\n");
    return 1;
  }

  bool match = XtalComp::compare(cell1, types1, pos1,
                                 cell2, types2, pos2,
                                 transform, cartTol, angleTol);

  printf("XtalComp Results\n");
  printf("Result:\n");
  printf("Using a cartesian tolerance of %f and an angular tolerance of %f...\n",
         cartTol, angleTol);
  printf("The structures %s match!\n", (match) ? "DO" : "do NOT");
  if (match) { // Print transform
    printf("Transformation matrix:\n");
    printf("|--%10s--%10s--%10s--%10s--|\n", "----------", "----------",
           "----------", "----------");
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[0*4+0], transform[0*4+1], transform[0*4+2], transform[0*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[1*4+0], transform[1*4+1], transform[1*4+2], transform[1*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[2*4+0], transform[2*4+1], transform[2*4+2], transform[2*4+3]);
    printf("|  %+10.5f  %+10.5f  %+10.5f  %+10.5f  |\n",
           transform[3*4+0], transform[3*4+1], transform[3*4+2], transform[3*4+3]);
    printf("|--%10s--%10s--%10s--%10s--|\n", "----------", "----------",
           "----------", "----------");
    }

  printf("Input structures:\n") ;
  PRINT_DIV;
  printf("| %-35s | %-35s |\n",
         "First cell matrix (row vectors)",
         "Second cell matrix (row vectors)");
  PRINT_DIV;
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[0][0], cell1[0][1], cell1[0][2], "",
         cell2[0][0], cell2[0][1], cell2[0][2], "");
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[1][0], cell1[1][1], cell1[1][2], "",
         cell2[1][0], cell2[1][1], cell2[1][2], "");
  printf("| %9.5f %9.5f %9.5f %5s | %9.5f %9.5f %9.5f %5s |\n",
         cell1[2][0], cell1[2][1], cell1[2][2], "",
         cell2[2][0], cell2[2][1], cell2[2][2], "");
  PRINT_DIV;
  printf("| %-35s | %-35s |\n",
         "type: fractional coordinate",
         "type: fractional coordinate");
  PRINT_DIV;
  for (int i = 0; i < types1.size(); ++i) {
    printf("| %3hu: %9.5f %9.5f %9.5f %0s | %3hu: %9.5f %9.5f %9.5f %0s |\n",
           types1[i], pos1[i][0], pos1[i][1], pos1[i][2], "",
           types2[i], pos2[i][0], pos2[i][1], pos2[i][2], "");
  }
  PRINT_DIV;

  return 0;
}




