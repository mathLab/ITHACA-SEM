<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit, Adiabatic wall</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_Adiabatic.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_Adiabatic.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-7">0.279194</value>
            <value variable="rhou" tolerance="1e-6">26.3093</value>
            <value variable="rhov" tolerance="1e-8">0.337764</value>
            <value variable="E" tolerance="1e-12">444.252</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.249122</value>
            <value variable="rhou" tolerance="1e-12">31.5414</value>
            <value variable="rhov" tolerance="1e-8">0.745257</value>
            <value variable="E" tolerance="1e-12">620.575</value>
        </metric>
    </metrics>
</test>
