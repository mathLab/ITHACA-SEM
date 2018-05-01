<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>New RiemannInvariant boundary to keep velocity at inflow fixed at u_inf</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>--use-scotch RiemannBoundary.xml</parameters>
    <files>
        <file description="Session File"> RiemannBoundary.xml</file>
        <file description="Restart File"> RiemannInvariantBC_Original.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12"> 0.989112</value>
            <value variable="rhou" tolerance="1e-12">366.648</value>
            <value variable="rhov" tolerance="1e-12">71.0173</value>
            <value variable="E" tolerance="1e-12">1.03847e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.320136</value>
            <value variable="rhou" tolerance="1e-12">152.635</value>
            <value variable="rhov" tolerance="1e-12">80.8759</value>
            <value variable="E" tolerance="1e-12">281783</value>
        </metric>
    </metrics>
</test>


