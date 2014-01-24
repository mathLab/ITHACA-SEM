<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Subsonic Cylinder, mixed bcs, FRDG advection and LFRDG diffusion, GLL_LAGRANGE</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">CylinderSubsonic_NS_FRDG_LFRDG_GLL_LAGRANGE_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">3.89457</value>
            <value variable="rhou" tolerance="1e-12">1425.78</value>
            <value variable="rhov" tolerance="1e-12">17.347</value>
            <value variable="rhow" tolerance="1e-12">17.3092</value>
            <value variable="E" tolerance="1e-12">4.43372e+06</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.29138</value>
            <value variable="rhou" tolerance="1e-12">86.3215</value>
            <value variable="rhov" tolerance="1e-12"> 14.0127</value>
            <value variable="rhow" tolerance="1e-12">1.00005</value>
            <value variable="E" tolerance="1e-12">275768</value>
        </metric>
    </metrics>
</test>


