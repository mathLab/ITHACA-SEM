<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P8.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P8.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">6.341e-05</value>
            <value variable="rhou" tolerance="1e-12">0.000156986</value>
            <value variable="rhov" tolerance="1e-12">0.0148192</value>
            <value variable="E" tolerance="1e-12">0.0082784</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00057786</value>
            <value variable="rhou" tolerance="1e-12">0.00181699</value>
            <value variable="rhov" tolerance="1e-12">0.00817004</value>
            <value variable="E" tolerance="1e-12">0.0093858</value>
        </metric>
    </metrics>
</test>


