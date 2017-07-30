<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, quasi 1D nozzle, stagnation inflow bc</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Nozzle_Quasi1D_P6.xml</parameters>
    <files>
        <file description="Session File">Nozzle_Quasi1D_P6.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.86343</value>
            <value variable="rhou" tolerance="1e-12">27.0141</value>
            <value variable="E" tolerance="1e-12">787959</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.26257</value>
            <value variable="rhou" tolerance="1e-12">50.4078</value>
            <value variable="E" tolerance="1e-12">260694</value>
        </metric>
    </metrics>
</test>


