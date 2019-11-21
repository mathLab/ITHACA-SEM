<?xml version="1.0" encoding="utf-8"?>
<test>
    <description> NekMesh with Star CCM input and STEP file to project CAD </description>
    <executable>NekMesh</executable>
    <parameters> -m projectcad:file=Ahmed_body.STEP -m jac:list Ahmed_body.ccm Ahmed_body.xml </parameters>
    <files>
        <file description="Input File">Ahmed_body.ccm</file>
        <file description="Input File">Ahmed_body.STEP</file>
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
