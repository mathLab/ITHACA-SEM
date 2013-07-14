<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and LDG diffusion, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_LDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_LDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.0885069</value>
            <value variable="rhou" tolerance="1e-12">60.3481</value>
            <value variable="rhov" tolerance="1e-8">0.215084</value>
            <value variable="E" tolerance="1e-12">4869.21</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0746837</value>
            <value variable="rhou" tolerance="1e-12">60.8017</value>
            <value variable="rhov" tolerance="1e-8">0.263229</value>
            <value variable="E" tolerance="1e-12">4404.81</value>
        </metric>
    </metrics>
</test>


