#include "guineapigCPP.h"
#include "particleBeamCPP.h"
#include "option_args.h"

int mainGuineapig(std::vector<std::string> input_arguments)
{
  // There is probably a nicer way to do this...
   std::string arg1 = input_arguments[0];
   std::string arg2 = input_arguments[1];
   std::string arg3 = input_arguments[2];

   GUINEA guinee((char*)(arg1.c_str()));
   guinee.run((char*)(arg2.c_str()),(char*)(arg3.c_str()));
   return 0;
}

           
  
int main (int argc,char *argv[])
{
  std::vector<std::string> input_arguments = read_args(argc,argv);
  return mainGuineapig(input_arguments);
}
