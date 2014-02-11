<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, FRSD advection and LFRSD diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Couette_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000391612</value>
            <value variable="rhou" tolerance="1e-12">48.1285</value>
            <value variable="rhov" tolerance="1e-12">0.143446</value>
            <value variable="rhow" tolerance="1e-12"> 6.98852e-06</value>
            <value variable="E" tolerance="1e-12">17519.5</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00164134</value>
            <value variable="rhou" tolerance="1e-12">83.3449</value>
            <value variable="rhov" tolerance="1e-12">0.5882</value>
            <value variable="rhow" tolerance="1e-12">4.62189e-05</value>
            <value variable="E" tolerance="1e-12">18878.3</value>
        </metric>
    </metrics>
</test>


