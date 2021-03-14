<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>NS, Couette flow, mixed bcs, WeakDG advection and IP diffusion, Implicit</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</parameters>
    <processes>4</processes>
    <files>
        <file description="Session File">Couette_WeakDG_IP_MODIFIED_IM_PRECON.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-4">3.40370e-04</value>
            <value variable="rhou" tolerance="2e-2">1.57590e-01</value>
            <value variable="rhov" tolerance="3e-4">3.47046e-02</value>
            <value variable="E" tolerance="2e-3">3.24624e+01</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="5e-4">2.78404e-04</value>
            <value variable="rhou" tolerance="2e-2">1.32515e-01</value>
            <value variable="rhov" tolerance="3e-4">4.35072e-02</value>
            <value variable="E" tolerance="5e-1">3.37731e+01</value>
        </metric>

    </metrics>
</test>
