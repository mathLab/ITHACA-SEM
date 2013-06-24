<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, FRSD, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_FRSD_MODIFIED_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_FRSD_MODIFIED_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">8.44263e-05</value>
            <value variable="rhou" tolerance="1e-12">0.000284056</value>
            <value variable="rhov" tolerance="1e-12">0.000274199</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-12">0.000804203</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12"> 1.10842e-06</value>
            <value variable="rhou" tolerance="1e-12"> 2.87238e-06</value>
            <value variable="rhov" tolerance="1e-12">4.24192e-06</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-12">1.25186e-05</value>
        </metric>
    </metrics>
</test>


