<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=8 GAUSS</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_P8_GAUSS.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_P8_GAUSS.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.40009e-07</value>
            <value variable="rhou" tolerance="1e-12">3.57402e-07</value>
            <value variable="rhov" tolerance="1e-12">3.78317e-07</value>
            <value variable="E" tolerance="1e-12">1.26238e-06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">3.77656e-07</value>
            <value variable="rhou" tolerance="1e-12">1.34087e-06</value>
            <value variable="rhov" tolerance="1e-12">1.39239e-06</value>
            <value variable="E" tolerance="1e-12">5.19969e-06</value>
        </metric>
    </metrics>
</test>


