#include "Process.h"

#include <set>
#include <iostream>

Process::Process(std::string _name , Color_t _color)
	: name(_name),
	  color(_color)
{
	for ( const auto& pol : {LL,LR,RL,RR} )
	{
		nSelected.insert( { pol , 0 } ) ;
		nTotal.insert( { pol , 0 } ) ;
	}
}

int Process::getProcess(const int processID, const bool isSignal)
{
	const std::set<int> signal = { 402007 , 402008 } ;
	const std::set<int> higgsstrahlung = { 402001 , 402002 , 402003 , 402004 , 402005 , 402006 , 402009 , 402010 , 402011 , 402012 } ;
	const std::set<int> _2f = { 500006 , 500008 , 500010 , 500012 } ;
	const std::set<int> _4f = { 500070 , 500072 , 500074 , 500076 , 500082 , 500084 , 500098 , 500100 , 500101 , 500102 , 500103 , 500104 , 500105 , 500106 , 500107 , 500108 , 500110 , 500112 , 					     500113 , 500114 , 500115 , 500116 , 500117 , 500118 , 500119 , 500120 , 500122 , 500124 } ;

	if (isSignal && !signal.count(processID))
		throw;
		
	if ( signal.count(processID) )
	{
		if (isSignal)
			return 0 ;
		else
			return 1;
	}
	else if ( higgsstrahlung.count(processID) )
		return 2 ;
	else if ( _2f.count(processID) )
		return 3 ;
	else if ( _4f.count(processID) )
		return 4 ;
	else
	{
		std::cout << "Process::Error in getProcess : processID = " << processID << std::endl ;
		throw ;
	}
}
