<?xml version="1.0" encoding="utf-8"?>
<test>    
    <description>Euler 1D P=3, WeakDG, MODIFIED</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>Euler1D.xml</parameters>
    <files>
        <file description="Session File">Euler1D.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">1.98838e-06</value>
            <value variable="rhou" tolerance="1e-12">0.00067684</value>
            <value variable="E" tolerance="1e-12">0.575708</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">1.98524e-05</value>
            <value variable="rhou" tolerance="1e-12">0.00675771</value>
            <value variable="E" tolerance="1e-12">5.74799</value>
        </metric>
    </metrics>
</test>

