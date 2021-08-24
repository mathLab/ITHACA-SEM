<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=11, periodic bcs, FRHU, GLL_LAGRANGE_SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex_FRHU_SEM.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex_FRHU_SEM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0457922</value>
            <value variable="rhou" tolerance="1e-12">0.0791701</value>
            <value variable="rhov" tolerance="1e-12">0.0732046</value>
            <value variable="E" tolerance="1e-12">0.181007</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0495651</value>
            <value variable="rhou" tolerance="1e-12">0.0949787</value>
            <value variable="rhov" tolerance="1e-12">0.136381</value>
            <value variable="E" tolerance="1e-12">0.269242</value>
        </metric>
    </metrics>
</test>


