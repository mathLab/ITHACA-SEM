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
            <value variable="rho" tolerance="1e-12">3.86521</value>
            <value variable="rhou" tolerance="1e-12">24.6372</value>
            <value variable="E" tolerance="1e-12">787766</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.2538</value>
            <value variable="rhou" tolerance="1e-12">45.0691</value>
            <value variable="E" tolerance="1e-12">260066</value>
        </metric>
    </metrics>
</test>