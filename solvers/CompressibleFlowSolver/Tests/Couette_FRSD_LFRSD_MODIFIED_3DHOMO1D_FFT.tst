<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRSD advection and LFRSD diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000391074</value>
            <value variable="rhou" tolerance="1e-12">48.1285</value>
            <value variable="rhov" tolerance="1e-12">0.143445</value>
            <value variable="rhow" tolerance="1e-12"> 6.96136e-06</value>
            <value variable="E" tolerance="1e-12">17519.5</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0016255</value>
            <value variable="rhou" tolerance="1e-12">83.3449</value>
            <value variable="rhov" tolerance="1e-12">0.588196</value>
            <value variable="rhow" tolerance="1e-12">4.16796e-05</value>
            <value variable="E" tolerance="1e-12">18878.4</value>
        </metric>
    </metrics>
</test>


