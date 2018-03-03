<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, axisymmetric pipe flow with low Mach number</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>PipeFlow_NSAxisym.xml</parameters>
    <files>
        <file description="Session File"> PipeFlow_NSAxisym.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-9">1.64331e-05</value>
            <value variable="rhou" tolerance="1e-12">0.00300872</value>
            <value variable="rhov" tolerance="1e-12">0.00266976</value>
            <value variable="E" tolerance="1e-12">4.484</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">2.00604e-05</value>
            <value variable="rhou" tolerance="1e-12">0.00264886</value>
            <value variable="rhov" tolerance="1e-12">0.0031863</value>
            <value variable="E" tolerance="1e-12">5.412</value>
        </metric>
    </metrics>
</test>

