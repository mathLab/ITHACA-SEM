<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Ringleb Flow P=3</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>RinglebFlow_P3.xml</parameters>
    <files>
        <file description="Session File">RinglebFlow_P3.xml</file>
        <file description="Restart File">RinglebFlow_P3.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000570231</value>
            <value variable="rhou" tolerance="1e-12">0.0011389</value>
            <value variable="rhov" tolerance="1e-12">0.000470153</value>
            <value variable="E" tolerance="1e-12">0.00100199</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00546575</value>
            <value variable="rhou" tolerance="1e-12">0.0171883</value>
            <value variable="rhov" tolerance="1e-12">0.00163612</value>
            <value variable="E" tolerance="1e-12">0.00892688</value>
        </metric>
    </metrics>
</test>


