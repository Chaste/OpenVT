#ifndef TEST05CELLSORTING_HPP_
#define TEST05CELLSORTING_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "PetscSetupAndFinalize.hpp"



class Test05CellSorting: public AbstractCellBasedWithTimingsTestSuite
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

#endif /* TEST05CELLSORTING_HPP_ */
