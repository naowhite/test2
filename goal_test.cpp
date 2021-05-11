#include <windows.h>
#include <iostream>
#include <list>

using namespace std;

enum class GoalNodeStatus
{
	Stopped,
	Running,
	Succeeded,
	Failed,
};

template <class EntityType, class Real=double>
struct GoalUpdateArgs
{
	GoalUpdateArgs(EntityType& entity, const Real& dt) : entity(entity), dt(dt) {}
	EntityType& entity;
	const Real& dt;
};

template <class EntityType, class Real=double>
class GoalNode
{
public:
	typedef GoalNode<EntityType, Real> NodeType;
	typedef GoalUpdateArgs<EntityType, Real> GoalUpdateArgsType;

	GoalNode() : _status(GoalNodeStatus::Stopped) {}
	virtual ~GoalNode() = default;
	virtual void OnStart(EntityType&) {}
	virtual void OnEnd(EntityType&) {}
	virtual GoalNodeStatus OnUpdate(GoalUpdateArgsType&) abstract;
	virtual void RemoveAllSubGoals(EntityType& entity) {}

	bool IsStopped() const { return _status == GoalNodeStatus::Stopped; }
	
	void Start(EntityType& entity) 
	{
		OnStart(entity);
		_status = GoalNodeStatus::Running;
	}

	void End(EntityType& entity)
	{
		OnEnd(entity);
		_status = GoalNodeStatus::Stopped;
	}

	GoalNodeStatus Update(GoalUpdateArgsType& args) 
	{
		_status = OnUpdate(args);
		return _status;
	}

protected:
	GoalNodeStatus _status;
};


template <class EntityType, class Real=double>
class GoalComposite : public GoalNode<EntityType,Real>
{
public:
	typedef GoalNode<EntityType, Real> NodeType;
	typedef GoalUpdateArgs<EntityType, Real> GoalUpdateArgsType;

	virtual ~GoalComposite() = default;

	template <class UNodeType>
	UNodeType* PushBackSubGoal(UNodeType* node)
	{
		_subgoals.push_back(node);
		return node;
	}

	template <class UNodeType>
	UNodeType* PushFrontSubGoal(UNodeType* node)
	{
		_subgoals.push_front(node);
		return node;
	}
	
private:
	void RemoveAllSubGoals(EntityType& entity)
	{
		for (NodeType* subgoal : _subgoals)
		{
			subgoal->RemoveAllSubGoals(entity);
			if (!subgoal->IsStopped())
			{
				subgoal->End(entity);
			}
		}
		_subgoals.clear();
	}

public:
	virtual GoalNodeStatus OnUpdate(GoalUpdateArgsType& args) override final
	{
		if (this->IsStopped())
		{
			this->Start(args.entity);
		}

		while (!_subgoals.empty())
		{
			NodeType& front = *_subgoals.front();
			if (front.IsStopped())
			{
				front.Start(args.entity);
			}

			GoalNodeStatus status = front.Update(args);
			switch (status)
			{
			case GoalNodeStatus::Succeeded:
				front.End(args.entity);
				_subgoals.pop_front();
				OnSubGoalSucceeded(args.entity);
				break;

			case GoalNodeStatus::Failed:
				RemoveAllSubGoals(args.entity);
				this->End(args.entity);
				OnSubGoalFailed(args.entity);
				return GoalNodeStatus::Failed;

			case GoalNodeStatus::Running:
				return GoalNodeStatus::Running;
			}
		}

		this->End(args.entity);
		OnAllSubGoalsSucceeded(args.entity);
		return GoalNodeStatus::Succeeded;
	}

	virtual void OnSubGoalSucceeded(EntityType&)
	{}

	virtual void OnSubGoalFailed(EntityType&)
	{
	}

	virtual void OnAllSubGoalsSucceeded(EntityType&)
	{
	}

private:
	std::list<NodeType*> _subgoals;
};

template <class EntityType, class Real=double>
class GoalTaskNull : public GoalNode<EntityType,Real>
{
public:
	typedef GoalUpdateArgs<EntityType, Real> GoalUpdateArgsType;
	virtual ~GoalTaskNull() {}
	virtual GoalNodeStatus OnUpdate(GoalUpdateArgsType&) { return GoalNodeStatus::Succeeded; }
};

// 以下、ユーザーコード

class TestEntity;

class GoalThink : public GoalComposite<TestEntity>
{
public:
	typedef GoalComposite<TestEntity> BaseType;
	typedef TestEntity EntityType;
	GoalThink();

	virtual void OnStart(EntityType& entity)
	{
		cout << "GoalThink::OnStart" << endl;
	}

	virtual void OnEnd(EntityType& entity)
	{
		cout << "GoalThink::OnEnd" << endl;
	}


	virtual void OnSubGoalSucceeded(EntityType& entity) override
	{
		cout << "GoalThink::OnSubGoalSucceeded" << endl;
	}

	virtual void OnSubGoalFailed(EntityType& entity) override;
	virtual void OnAllSubGoalsSucceeded(EntityType& entity) override;
};

class GoalThinkSub : public GoalComposite<TestEntity>
{
public:
	typedef TestEntity EntityType;
	GoalThinkSub();

	virtual void OnSubGoalSucceeded(EntityType& entity) override
	{
		cout << "GoalThinkSub::OnSubGoalSucceeded" << endl;
	}
};


class GoalTaskCreatePath : public GoalNode<TestEntity>
{
public:
	typedef TestEntity EntityType;

	GoalTaskCreatePath() {}

	virtual void OnStart(EntityType& entity) override
	{
		cout << "GoalTaskCreatePath::OnStart" << endl;
	}

	virtual void OnEnd(EntityType& entity) override
	{
		cout << "GoalTaskCreatePath::OnEnd" << endl;
	}

	virtual GoalNodeStatus OnUpdate(GoalUpdateArgsType& args) override
	{
		cout << "GoalTaskCreatePath::OnUpdate" << endl;
		return GoalNodeStatus::Succeeded;
		//return GoalNodeStatus::Failed;
	}
};

class GoalTaskFollowPath : public GoalNode<TestEntity>
{
public:
	typedef TestEntity EntityType;

	GoalTaskFollowPath() {}

	virtual void OnStart(EntityType& entity) override
	{
		cout << "GoalTaskFollowPath::OnStart" << endl;
	}

	virtual void OnEnd(EntityType& entity) override
	{
		cout << "GoalTaskFollowPath::OnEnd" << endl;
	}

	virtual GoalNodeStatus OnUpdate(GoalUpdateArgsType& args) override
	{
		cout << "GoalTaskFollowPath::OnUpdate" << endl;
		return GoalNodeStatus::Succeeded;
		//return GoalNodeStatus::Failed;
	}
};

class TestEntity
{
public:
	typedef GoalUpdateArgs<TestEntity> GoalUpdateArgsType;

	TestEntity();

	void Update()
	{
		GoalUpdateArgsType args(*this, 0.1);
		mGoalThink.Update(args);
	}

	void SetGoal();

	GoalThink mGoalThink;
	GoalThinkSub mGoalThinkSub;
	GoalTaskCreatePath mGoalTaskCreatePath;
	GoalTaskFollowPath mGoalTaskFollowPath;
	GoalTaskCreatePath mGoalTaskCreatePath2;
	GoalTaskFollowPath mGoalTaskFollowPath2;

};

TestEntity::TestEntity()
{
	SetGoal();
}

void TestEntity::SetGoal()
{
	mGoalThink.PushBackSubGoal(&mGoalTaskCreatePath);
	mGoalThink.PushBackSubGoal(&mGoalTaskFollowPath);
	auto pGoalThinkSub = mGoalThink.PushBackSubGoal(&mGoalThinkSub);
	pGoalThinkSub->PushBackSubGoal(&mGoalTaskCreatePath2);
	pGoalThinkSub->PushBackSubGoal(&mGoalTaskFollowPath2);
}

GoalThink::GoalThink()
{
}

GoalThinkSub::GoalThinkSub()
{
}


void GoalThink::OnSubGoalFailed(EntityType& entity)
{
	cout << "GoalThink::OnSubGoalFailed" << endl;
	entity.mGoalThink.PushBackSubGoal(&entity.mGoalTaskCreatePath);
	entity.mGoalThink.PushBackSubGoal(&entity.mGoalTaskFollowPath);
}

void GoalThink::OnAllSubGoalsSucceeded(EntityType& entity)
{
	cout << "GoalThink::OnAllSubGoalsSucceeded" << endl;
	entity.SetGoal();
}


int main()
{
	TestEntity a;

	while (true)
	{
		cout << "-- a.Update() --" << endl;
		a.Update();

		::Sleep(500);
	}

	return 0;
}