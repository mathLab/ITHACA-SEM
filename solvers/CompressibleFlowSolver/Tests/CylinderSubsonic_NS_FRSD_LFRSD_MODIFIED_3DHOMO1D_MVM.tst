<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, FRSD advection and LFRSD diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_FRSD_LFRSD_MODIFIED_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.89459</value>
            <value variable="rhou" tolerance="1e-12">1425.75</value>
            <value variable="rhov" tolerance="1e-12">17.4952</value>
            <value variable="rhow" tolerance="1e-12">17.3092</value>
            <value variable="E" tolerance="1e-12">4.43372e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.365224</value>
            <value variable="rhou" tolerance="1e-12">87.9294</value>
            <value variable="rhov" tolerance="1e-8">28.339</value>
            <value variable="rhow" tolerance="1e-12">1.00018</value>
            <value variable="E" tolerance="1e-12">297764</value>
        </metric>
    </metrics>
</test>


