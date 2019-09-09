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
            <value variable="rho" tolerance="1e-12">0.000570156</value>
            <value variable="rhou" tolerance="1e-12">0.00113889</value>
            <value variable="rhov" tolerance="1e-12">0.000470214</value>
            <value variable="E" tolerance="1e-12">0.00100227</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.00546578</value>
            <value variable="rhou" tolerance="1e-12">0.0171878</value>
            <value variable="rhov" tolerance="1e-12">0.00163611</value>
            <value variable="E" tolerance="1e-12">0.008927</value>
        </metric>
    </metrics>
</test>


