<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Linearized Channel Flow P=11</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChannelSpongeLNSE.xml</parameters>
    <files>
        <file description="Session File">ChannelSpongeLNSE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">0.000360221</value>
            <value variable="v" tolerance="1e-12">0.000186555</value>
	    <value variable="p" tolerance="1e-12">0.000232531</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">0.000924622</value>
            <value variable="v" tolerance="1e-12">0.00057635</value>
	    <value variable="p" tolerance="1e-12">0.000671219</value>
        </metric>
    </metrics>
</test>


