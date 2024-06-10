#ifndef TEST03WOUNDHEALING_HPP_
#define TEST03WOUNDHEALING_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "PetscSetupAndFinalize.hpp"



class Test03WoundHealing: public AbstractCellBasedWithTimingsTestSuite
{
private:

public:

    /*
     */
   void TestNotMuch()
   {
        TS_ASSERT_THROWS_NOWT();
    }
};

#endif /* TEST03WOUNDHEALING_HPP_ */
