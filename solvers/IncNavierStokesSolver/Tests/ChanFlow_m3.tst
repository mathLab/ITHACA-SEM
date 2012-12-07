<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=3</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChanFlow_m3.xml</parameters>
    <files>
        <file description="Session File">ChanFlow_m3.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.49822e-16</value>
            <value variable="v" tolerance="1e-12">0</value>
	    <value variable="p" tolerance="1e-12">8.05332e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">3.60822e-15</value>
            <value variable="v" tolerance="1e-12">2.70499e-16</value>
	    <value variable="p" tolerance="1e-12">6.50591e-14</value>
        </metric>
    </metrics>
</test>


