#ifndef TEST01PERSISTENTRANDOMWALK_HPP_
#define TEST01PERSISTENTRANDOMWALK_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "PetscSetupAndFinalize.hpp"



class Test01PersistentRandomWalk: public AbstractCellBasedWithTimingsTestSuite
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

#endif /* TEST01PERSISTENTRANDOMWALK_HPP_ */
