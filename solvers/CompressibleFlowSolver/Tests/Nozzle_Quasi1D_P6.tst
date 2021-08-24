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
            <value variable="rho" tolerance="1e-12">3.87376</value>
            <value variable="rhou" tolerance="1e-12">3.89851</value>
            <value variable="E" tolerance="1e-12">790564</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.22627</value>
            <value variable="rhou" tolerance="1e-12">2.95978</value>
            <value variable="E" tolerance="1e-12">250372</value>
        </metric>
    </metrics>
</test>