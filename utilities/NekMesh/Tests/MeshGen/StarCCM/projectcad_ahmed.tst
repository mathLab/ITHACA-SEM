<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh with Star CCM input and STEP file to project CAD </description>
    <executable>NekMesh</executable>
    <parameters> -m projectcad:file=projectcad_ahmed.stp -m jac:list projectcad_ahmed.ccm projectcad_ahmed-out.xml:xml:test </parameters>
    <files>
        <file description="Input File">projectcad_ahmed.ccm</file>
        <file description="Input File">projectcad_ahmed.stp</file>
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
