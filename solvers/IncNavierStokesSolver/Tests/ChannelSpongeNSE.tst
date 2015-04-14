<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Channel Flow P=11</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>ChannelSpongeNSE.xml</parameters>
    <files>
        <file description="Session File">ChannelSpongeNSE.xml</file>
    </files>
    <metrics>
        <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">2.589</value>
            <value variable="v" tolerance="1e-12">0.000171334</value>
	    <value variable="p" tolerance="1e-12">0.0128673</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">1</value>
            <value variable="v" tolerance="1e-12">0.000630178</value>
            <value variable="p" tolerance="1e-12">0.00746041</value>
        </metric>
    </metrics>
</test>


