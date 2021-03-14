<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-5">3.35211e-04</value>
            <value variable="rhou" tolerance="2e-3">1.56391e-01</value>
            <value variable="rhov" tolerance="2e-5">3.46942e-02</value>
            <value variable="E" tolerance="5e-4">3.24621e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="5e-5">3.00847e-04</value>
            <value variable="rhou" tolerance="2e-3">1.31175e-01</value>
            <value variable="rhov" tolerance="5e-6">4.35033e-02</value>
            <value variable="E" tolerance="5e-2">3.37993e+01</value>
        </metric>

    </metrics>
</test>
