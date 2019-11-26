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
            <value variable="rho" tolerance="1e-9">1.71551e-05</value>
            <value variable="rhou" tolerance="1e-8">0.00300236</value>
            <value variable="rhov" tolerance="1e-8">0.00267845</value>
            <value variable="E" tolerance="1e-4">4.6884</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-9">2.08147e-05</value>
            <value variable="rhou" tolerance="1e-8">0.002665</value>
            <value variable="rhov" tolerance="1e-8">0.00312463</value>
            <value variable="E" tolerance="1e-5">5.62463</value>
        </metric>
    </metrics>
</test>

