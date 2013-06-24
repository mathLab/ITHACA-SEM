<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRDG advection and LFRDG diffusion, GLL_LAGRANGE</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_FFT.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000483463</value>
            <value variable="rhou" tolerance="1e-12">48.1265</value>
            <value variable="rhov" tolerance="1e-12">0.177365</value>
            <value variable="rhow" tolerance="1e-12">7.17934e-06</value>
            <value variable="E" tolerance="1e-12">17518.5</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00207096</value>
            <value variable="rhou" tolerance="1e-12">83.3318</value>
            <value variable="rhov" tolerance="1e-12">0.752386</value>
            <value variable="rhow" tolerance="1e-12">4.41936e-05</value>
            <value variable="E" tolerance="1e-12">18722.8</value>
        </metric>
    </metrics>
</test>


