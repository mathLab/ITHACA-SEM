<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Kovasznay Flow P=8</description>
    <executable>IncNavierStokesSolver</executable>
    <parameters>KovaFlow_m8.xml</parameters>
    <files>
        <file description="Session File">KovaFlow_m8.xml</file>
	<file description="Session File">KovaFlow_m8.rst</file>
	<file description="Session File">KovaFlow_m8_avg.rst</file>
    </files>
    <metrics>
      <metric type="L2" id="1">
            <value variable="u" tolerance="1e-12">4.59041e-05</value>
            <value variable="v" tolerance="1e-12">0.000172901</value>
	    <value variable="p" tolerance="1e-12">4.0608e-05</value>
        </metric>
        <metric type="Linf" id="2">
            <value variable="u" tolerance="1e-12">7.46824e-05</value>
            <value variable="v" tolerance="1e-12">0.000224474</value>
	    <value variable="p" tolerance="1e-12">0.000614459</value>
        </metric>
    </metrics>
</test>
