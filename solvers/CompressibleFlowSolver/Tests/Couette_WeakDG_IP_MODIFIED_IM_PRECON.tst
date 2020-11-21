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
            <value variable="rho" tolerance="1e-8">0.00034037</value>
            <value variable="rhou" tolerance="1e-6">0.15759</value>
            <value variable="rhov" tolerance="1e-7">0.0347046</value>
            <value variable="E" tolerance="1e-4">32.4624</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-9">0.000278404</value>
            <value variable="rhou" tolerance="1e-6">0.132515</value>
            <value variable="rhov" tolerance="1e-7">0.0435072</value>
            <value variable="E" tolerance="1e-4">33.7731</value>
        </metric>

    </metrics>
</test>
