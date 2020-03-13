<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler Isentropic Vortex P=4, FRDG, GLL_LAGRANGE</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>IsentropicVortex16_FRDG_GLL_LAGRANGE_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">IsentropicVortex16_FRDG_GLL_LAGRANGE_3DHOMO1D_FFT.xml</file>
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
            <value variable="rho" tolerance="1e-12">0.000212639</value>
            <value variable="rhou" tolerance="1e-12"> 0.000904271</value>
            <value variable="rhov" tolerance="1e-12"> 0.000819096</value>
            <value variable="rhow" tolerance="1e-12">0</value>
            <value variable="E" tolerance="1e-12">0.00240946</value>
        </metric>
    </metrics>
</test>


