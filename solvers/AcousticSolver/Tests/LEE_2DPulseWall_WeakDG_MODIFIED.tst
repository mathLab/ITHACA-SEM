<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>desc P=400</description>
    <executable>AcousticSolver</executable>
    <parameters>LEE_2DPulseWall_WeakDG_MODIFIED.xml</parameters>
    <files>
        <file description="Session File">LEE_2DPulseWall_WeakDG_MODIFIED.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="p" tolerance="1e-4">  6.7577</value>
            <value variable="rho" tolerance="1e-12">9.49956e-05</value>
            <value variable="rhou" tolerance="1e-7"> 0.0175868</value>
            <value variable="rhov" tolerance="1e-7"> 0.00688346</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="p" tolerance="1e-4">  13.659</value>
            <value variable="rho" tolerance="1e-12">0.000327788</value>
            <value variable="rhou" tolerance="1e-7"> 0.0338246</value>
            <value variable="rhov" tolerance="1e-7"> 0.0131226</value>
        </metric>
    </metrics>
</test>
