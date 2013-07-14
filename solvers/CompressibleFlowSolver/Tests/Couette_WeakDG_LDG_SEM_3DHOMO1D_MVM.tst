<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, SEM</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_SEM_3DHOMO1D_MVM.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_SEM_3DHOMO1D_MVM.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000397744</value>
            <value variable="rhou" tolerance="1e-12">48.1307</value>
            <value variable="rhov" tolerance="1e-12">0.14585</value>
             <value variable="rhow" tolerance="1e-12">9.0851e-06</value>
            <value variable="E" tolerance="1e-12">17520.4</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00140076</value>
            <value variable="rhou" tolerance="1e-12">83.3598</value>
            <value variable="rhov" tolerance="1e-12">0.505652</value>
             <value variable="rhow" tolerance="1e-12">2.87012e-05</value>
            <value variable="E" tolerance="1e-12">18957.4</value>
        </metric>
    </metrics>
</test>


