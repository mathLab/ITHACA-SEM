<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P3.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.00169229</value>
            <value variable="rhou" tolerance="1e-12">0.00308706</value>
            <value variable="rhov" tolerance="1e-12">0.00481169</value>
            <value variable="E" tolerance="1e-12">0.0107409</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00312134</value>
            <value variable="rhou" tolerance="1e-12">0.00603793</value>
            <value variable="rhov" tolerance="1e-12">0.0110506</value>
            <value variable="E" tolerance="1e-12">0.0268249</value>
        </metric>
    </metrics>
</test>


