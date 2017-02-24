<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, FRHU advection and LFRHU diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_FRHU_LFRHU_SEM_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.89457</value>
            <value variable="rhou" tolerance="1e-12">1425.78</value>
            <value variable="rhov" tolerance="1e-12">17.3413</value>
            <value variable="rhow" tolerance="1e-12">17.3092</value>
            <value variable="E" tolerance="1e-12">4.43372e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12"> 0.275036</value>
            <value variable="rhou" tolerance="1e-12">82.5652</value>
            <value variable="rhov" tolerance="1e-8">10.8489</value>
            <value variable="rhow" tolerance="1e-12">1.00004</value>
            <value variable="E" tolerance="1e-12">270824</value>
        </metric>
    </metrics>
</test>


