<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRDG advection and LFRDG diffusion, GLL_LAGRANGE</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Couette_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000483154</value>
            <value variable="rhou" tolerance="1e-12">48.1265</value>
            <value variable="rhov" tolerance="1e-12">0.177367</value>
            <value variable="rhow" tolerance="1e-12"> 7.87069e-06</value>
            <value variable="E" tolerance="1e-12">17518.9</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00207196</value>
            <value variable="rhou" tolerance="1e-12">83.332</value>
            <value variable="rhov" tolerance="1e-12">0.752414</value>
            <value variable="rhow" tolerance="1e-12">5.89513e-05</value>
            <value variable="E" tolerance="1e-12">18726.3</value>
        </metric>
    </metrics>
</test>


