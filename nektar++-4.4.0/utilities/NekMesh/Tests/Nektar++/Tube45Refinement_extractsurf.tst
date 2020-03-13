<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>Test BL refinement from XML file</description>
    <executable>NekMesh</executable>
    <parameters>-m jac:list -m extract:surf=2:detectbnd Tube45Refinement.xml Tube45Refined.xml:xml:test</parameters>
    <files>
        <file description="Input File">Tube45Refinement.xml</file>
    </files>
    <metrics>
        <metric type="regex" id="1">
            <regex>^Total negative Jacobians: (\d+)</regex>
            <matches>
                <match>
                    <field id="0">0</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>
