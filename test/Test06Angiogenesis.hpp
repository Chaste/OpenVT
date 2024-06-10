#ifndef TEST06ANGIOGENESIS_HPP_
#define TEST06ANGIOGENESIS_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "PetscSetupAndFinalize.hpp"



class Test06Angiogenesis: public AbstractCellBasedWithTimingsTestSuite
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

#endif /* TEST06ANGIOGENESIS_HPP_ */
