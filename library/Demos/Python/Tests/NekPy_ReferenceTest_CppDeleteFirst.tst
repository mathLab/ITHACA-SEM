<?xml version="1.0" encoding="utf-8" ?>
<test>
	<description> Memory test for NekPy: Python to C++ array conversion </description>
	<executable python="true"> NekPy_ReferenceTest_CppDeleteFirst.py </executable>
	<parameters></parameters>
	<files>
        <file description="Session File">newsquare_2x2.xml</file>        
    </files>
	<metrics>
		<metric type="regex" id="1">
			<regex>^.*Test (.*)</regex>
            <matches>
                <match>
                    <field id="0">successful!</field>
                </match>
            </matches>			
		</metric>
	</metrics>
</test>