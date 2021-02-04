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
            <value variable="rho" tolerance="1e-7">0.0727295</value>
            <value variable="rhou" tolerance="1e-6">0.759621</value>
            <value variable="rhov" tolerance="1e-8">0.000927518</value>
            <value variable="E" tolerance="1e-12">7.68883</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0901262</value>
            <value variable="rhou" tolerance="1e-12">0.735283</value>
            <value variable="rhov" tolerance="1e-8">0.000934695</value>
            <value variable="E" tolerance="1e-12">7.97098</value>
        </metric>
    </metrics>
</test>


