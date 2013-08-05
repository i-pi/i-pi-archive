from common import TestSimulation

def test_lj_gas():
    ts = TestSimulation(input="../../test/lj/gas/input.xml", driver="../../drivers/driver.x")
    ts.run()
    # Test properties (e.g. latest positions/temperature etc)
