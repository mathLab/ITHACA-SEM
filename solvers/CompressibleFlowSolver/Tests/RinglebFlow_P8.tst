<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Euler, Ringleb Flow P=8</description>
    <executable>CompressibleFlowSolver</executable>
    <parameters>RinglebFlow_P8.xml</parameters>
    <files>
        <file description="Session File">RinglebFlow_P8.xml</file>
        <file description="Restart File">RinglebFlow_P8.rst</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="rho" tolerance="1e-12">0.000990531</value>
            <value variable="rhou" tolerance="1e-12">0.00326406</value>
            <value variable="rhov" tolerance="1e-12">0.000277794</value>
            <value variable="E" tolerance="1e-12">0.00165497</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="rho" tolerance="1e-12">0.0204517</value>
            <value variable="rhou" tolerance="1e-12">0.0407679</value>
            <value variable="rhov" tolerance="1e-12">0.00523929</value>
            <value variable="E" tolerance="1e-12">0.0323764</value>
        </metric>
    </metrics>
</test>


